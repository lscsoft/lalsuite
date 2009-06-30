/*
*  Copyright (C) 2007 Duncan Brown, Gareth Jones, Jolien Creighton
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
 * File Name: FindChirpBCVSpin.h
 *
 * Author: Brown, D. A. and Jones, G.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVSpinHV">
Author: Brown, D. A. and Jones, G
$Id$
</lalVerbatim>

<lalLaTeX>
\section{Header \texttt{FindChirpBCVSpin.h}}
\label{s:FindChirpBCVSpin.h}

Provides structures and functions to condition interferometer data
and generate binary inspiral chirps using the spinning BCV detection
template family.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpBCVSpin.h>
\end{verbatim}

\input{FindChirpBCVSpinHDoc}
</lalLaTeX>
#endif


#ifndef _FINDCHIRPBCVSPINH_H
#define _FINDCHIRPBCVSPINH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPBCVSPINH, "$Id$");

#if 0
<lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPBCVSPINH_ENULL 1
#define FINDCHIRPBCVSPINH_ENNUL 2
#define FINDCHIRPBCVSPINH_EALOC 3
#define FINDCHIRPBCVSPINH_ENUMZ 4
#define FINDCHIRPBCVSPINH_ESEGZ 5
#define FINDCHIRPBCVSPINH_EMISM 6
#define FINDCHIRPBCVSPINH_EDELT 7
#define FINDCHIRPBCVSPINH_EFLOW 8
#define FINDCHIRPBCVSPINH_EDYNR 9
#define FINDCHIRPBCVSPINH_EISTN 10
#define FINDCHIRPBCVSPINH_EDIVZ 11
#define FINDCHIRPBCVSPINH_EMAPX 12
#define FINDCHIRPBCVSPINH_EUAPX 13
#define FINDCHIRPBCVSPINH_ECLUW 14
#define FINDCHIRPBCVSPINH_MSGENULL "Null pointer"
#define FINDCHIRPBCVSPINH_MSGENNUL "Non-null pointer"
#define FINDCHIRPBCVSPINH_MSGEALOC "Memory allocation error"
#define FINDCHIRPBCVSPINH_MSGENUMZ "Invalid number of segments"
#define FINDCHIRPBCVSPINH_MSGESEGZ "Invalid number of points in segments"
#define FINDCHIRPBCVSPINH_MSGEMISM "Mismatch between number of points in segments"
#define FINDCHIRPBCVSPINH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPBCVSPINH_MSGEFLOW "Low frequency cutoff is negative"
#define FINDCHIRPBCVSPINH_MSGEDYNR "Dynamic range scaling is zero or negative"
#define FINDCHIRPBCVSPINH_MSGEISTN "Truncation of inverse power spectrum is negative"
#define FINDCHIRPBCVSPINH_MSGEDIVZ "Attempting to divide by zero"
#define FINDCHIRPBCVSPINH_MSGEMAPX "Mismatch in waveform approximant"
#define FINDCHIRPBCVSPINH_MSGEUAPX "Unknown approximant"
#define FINDCHIRPBCVSPINH_MSGECLUW "Unacceptable max-over-chirp clustering method for BCVSpin"

/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Types}

None.
</lalLaTeX>
#endif


#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpBCVSpinHV}}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\newpage\input{FindChirpBCVSpinDataC}
</lalLaTeX>
#endif

void
LALFindChirpBCVSpinData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpBCVSpinTemplateC}
</lalLaTeX>
#endif

void
LALFindChirpBCVSpinTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *thistmplt,
    FindChirpTmpltParams       *params,
    FindChirpDataParams        *fcDataParams
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpBCVSpinFilterC}
</lalLaTeX>
#endif

void
LALFindChirpBCVSpinFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    FindChirpDataParams        *fcDataParams
  );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPSPH_H */
