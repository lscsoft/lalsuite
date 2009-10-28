/*
*  Copyright (C) 2007 Darren Woods, Duncan Brown, Eirini Messaritaki, Gareth Jones, Jolien Creighton
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
 * File Name: FindChirpSP.h
 *
 * Author: Brown, D. A.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpSPHV">
Author: Brown, D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\section{Header \texttt{FindChirpSP.h}}
\label{s:FindChirpSP.h}

Provides structures and functions to condition interferometer data
and generate binary inspiral chirps using the stationary phase approximation.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpSP.h>
\end{verbatim}

\input{FindChirpSPHDoc}
</lalLaTeX>
#endif


#ifndef _FINDCHIRPSPH_H
#define _FINDCHIRPSPH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPSPH, "$Id$");

#if 0
<lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPSPH_ENULL 1
#define FINDCHIRPSPH_ENNUL 2
#define FINDCHIRPSPH_EALOC 3
#define FINDCHIRPSPH_ENUMZ 4
#define FINDCHIRPSPH_ESEGZ 5
#define FINDCHIRPSPH_EMISM 6
#define FINDCHIRPSPH_EDELT 7
#define FINDCHIRPSPH_EFLOW 8
#define FINDCHIRPSPH_EDYNR 9
#define FINDCHIRPSPH_EISTN 10
#define FINDCHIRPSPH_EDIVZ 11
#define FINDCHIRPSPH_EMAPX 12
#define FINDCHIRPSPH_EUAPX 13
#define FINDCHIRPSPH_EORDR 14
#define FINDCHIRPSPH_MSGENULL "Null pointer"
#define FINDCHIRPSPH_MSGENNUL "Non-null pointer"
#define FINDCHIRPSPH_MSGEALOC "Memory allocation error"
#define FINDCHIRPSPH_MSGENUMZ "Invalid number of segments"
#define FINDCHIRPSPH_MSGESEGZ "Invalid number of points in segments"
#define FINDCHIRPSPH_MSGEMISM "Mismatch between number of points in segments"
#define FINDCHIRPSPH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPSPH_MSGEFLOW "Low frequency cutoff is negative"
#define FINDCHIRPSPH_MSGEDYNR "Dynamic range scaling is zero or negative"
#define FINDCHIRPSPH_MSGEISTN "Truncation of inverse power spectrum is negative"
#define FINDCHIRPSPH_MSGEDIVZ "Attempting to divide by zero"
#define FINDCHIRPSPH_MSGEMAPX "Mismatch in waveform approximant"
#define FINDCHIRPSPH_MSGEUAPX "Unknown approximant"
#define FINDCHIRPSPH_MSGEORDR "Invalid post-Newtonian order"
/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Types}

None.
</lalLaTeX>
#endif


#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpSPHV}}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\newpage\input{FindChirpSPDataC}
</lalLaTeX>
#endif

void
LALFindChirpSPData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpSPTemplateC}
</lalLaTeX>
#endif

void
LALFindChirpSPTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPSPH_H */
