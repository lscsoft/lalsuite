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

