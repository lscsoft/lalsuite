/*
*  Copyright (C) 2007 Duncan Brown
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpBCVChisq.c
 *
 * Author: Anderson, W. G., and Brown, D. A., Spin-BCV-Modifications: Jones, G
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVSpinChisqCV">
Author: Anderson, W. G., and Brown D. A., Spin-BCV-Modifications: Jones, G
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpBCVSpinChisq.c}}
\label{ss:FindChirpBCVSpinChisq.c}

Module to implement the $\chi^2$ veto for the spinning BCV templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
% add this back in once the function is prototyped
%\input{FindChirpBCVSpinChisqCP}
\idx{LALFindChirpBCVSpinChisqVeto()}

\subsubsection*{Description}

The function \texttt{LALFindChirpBCVSpinChisqVeto()} perfoms a $\chi^2$ veto
on an entire data segment using the corresponding algorithm for the spinning
BCV templates, described below. On exit the vector \texttt{chisqVec} contains
the value $\chi^2(t_j)$ for the data segment.

\subsubsection*{Algorithm}

chisq algorithm here

\subsubsection*{Uses}
\begin{verbatim}
LALCreateReverseComplexFFTPlan()
LALDestroyComplexFFTPlan()
LALCCreateVector()
LALCDestroyVector()
LALCOMPLEX8VectorFFT()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpBCVSpinChisqCV}}
</lalLaTeX>
#endif

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include <lal/FindChirpBCVSpin.h>

NRCSID (FINDCHIRPBCVSPINCHISQC, "$Id$");

