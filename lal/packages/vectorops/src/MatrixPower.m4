/******** <lalVerbatim file="MatrixPowerCV"> ********
Author: Tibbits, M. M.
$Id$
********* </lalVerbatim> ********/

/******************************* <lalLaTeX file="MatrixPowerC">
\subsection{Module \texttt{MatrixPower.c}}
\label{s:MatrixPower.c}
This file is dedicated to reproducing the matlab function " \verb+.^+ ".  This file
has several declarations of the same function taking all forms of available
input.  This being said, I have yet to script the complex actions and their
counterparts.


\subsection*{Description}

This file is to help make the conversion from Matlab to c much earier.
In this file, we have created all of the versions of \verb+.^+ that we plan on
using.

\subsection*{Algorithms}

The algorithm is the same as it is in matlab.  The dot in front of an operator
in matlab signifies that if either or both of the operands are vectors, then
the operation will be carried out member by member.  For instance

vector a[25];
vector b[25];
vector c[25];

c = a \verb+.^+ b;

The result of this is:

\(c[0] =	a[0] ^{b[0]}\);
\(c[1] =	a[1] ^{b[1]}\);
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

NRCSID( MATLABMATRIXPOWC, "$Id$");

define(`TYPECODE',`D')
define(`TYPECODE2',`D')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`D')
define(`TYPECODE2',`S')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`D')
define(`TYPECODE2',`I2')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`D')
define(`TYPECODE2',`I4')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`D')
define(`TYPECODE2',`I8')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`S')
define(`TYPECODE2',`S')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`S')
define(`TYPECODE2',`I2')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`S')
define(`TYPECODE2',`I4')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`S')
define(`TYPECODE2',`I8')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`I2')
define(`TYPECODE2',`I2')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`I4')
define(`TYPECODE2',`I2')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`I4')
define(`TYPECODE2',`I4')
include(`MatrixPowerBase.m4')


define(`TYPECODE',`I8')
define(`TYPECODE2',`I2')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`I8')
define(`TYPECODE2',`I4')
include(`MatrixPowerBase.m4')

define(`TYPECODE',`I8')
define(`TYPECODE2',`I8')
include(`MatrixPowerBase.m4')
