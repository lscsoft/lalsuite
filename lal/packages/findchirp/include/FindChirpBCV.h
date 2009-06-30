/*
*  Copyright (C) 2007 Duncan Brown, Eirini Messaritaki, Jolien Creighton, Thomas Cokelaer
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
 * File Name: FindChirpBCV.h
 *
 * Author: Brown, D. A. and Messaritaki, E.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVHV">
Author: Brown, D. A. and Messaritaki, E.
$Id$
</lalVerbatim>

<lalLaTeX>
\section{Header \texttt{FindChirpBCV.h}}
\label{s:FindChirpBCV.h}

Provides structures and functions to condition interferometer data
and generate binary inspiral chirps using the BCV detection template
family.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpBCV.h>
\end{verbatim}

\input{FindChirpBCVHDoc}
</lalLaTeX>
#endif


#ifndef _FINDCHIRPBCVH_H
#define _FINDCHIRPBCVH_H

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


NRCSID (FINDCHIRPBCVH, "$Id$");

#if 0
<lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPBCVH_ENULL 1
#define FINDCHIRPBCVH_ENNUL 2
#define FINDCHIRPBCVH_EALOC 3
#define FINDCHIRPBCVH_ENUMZ 4
#define FINDCHIRPBCVH_ESEGZ 5
#define FINDCHIRPBCVH_EMISM 6
#define FINDCHIRPBCVH_EDELT 7
#define FINDCHIRPBCVH_EFLOW 8
#define FINDCHIRPBCVH_EDYNR 9
#define FINDCHIRPBCVH_EISTN 10
#define FINDCHIRPBCVH_EDIVZ 11
#define FINDCHIRPBCVH_EMAPX 12
#define FINDCHIRPBCVH_EUAPX 13
#define FINDCHIRPBCVH_EZNRM 14
#define FINDCHIRPBCVH_EQLEN 15
#define FINDCHIRPBCVH_ECLUW 16
#define FINDCHIRPBCVH_MSGENULL "Null pointer"
#define FINDCHIRPBCVH_MSGENNUL "Non-null pointer"
#define FINDCHIRPBCVH_MSGEALOC "Memory allocation error"
#define FINDCHIRPBCVH_MSGENUMZ "Invalid number of segments"
#define FINDCHIRPBCVH_MSGESEGZ "Invalid number of points in segments"
#define FINDCHIRPBCVH_MSGEMISM "Mismatch between number of points in segments"
#define FINDCHIRPBCVH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPBCVH_MSGEFLOW "Low frequency cutoff is negative"
#define FINDCHIRPBCVH_MSGEDYNR "Dynamic range scaling is zero or negative"
#define FINDCHIRPBCVH_MSGEISTN "Truncation of inverse power spectrum is negative"
#define FINDCHIRPBCVH_MSGEDIVZ "Attempting to divide by zero"
#define FINDCHIRPBCVH_MSGEMAPX "Mismatch in waveform approximant"
#define FINDCHIRPBCVH_MSGEUAPX "Unknown approximant"
#define FINDCHIRPBCVH_MSGEZNRM "No non-zero value assigned to one of a1, b1, b2"
#define FINDCHIRPBCVH_MSGEQLEN "params->qVec->length not equal to params->qVecBCV->length"
#define FINDCHIRPBCVH_MSGECLUW "Unacceptable max-over-chirp clustering method for BCV"
/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Types}

None.
</lalLaTeX>
#endif


#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpBCVHV}}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\newpage\input{FindChirpBCVDataC}
</lalLaTeX>
#endif

void
LALFindChirpBCVData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpBCVTemplateC}
</lalLaTeX>
#endif

void
LALFindChirpBCVTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *theTmplt,
    FindChirpTmpltParams       *params
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpBCVChisqC}
</lalLaTeX>
#endif

void
LALFindChirpBCVChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqInput        *inputBCV,
    FindChirpChisqParams       *params
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpBCVFilterC}
</lalLaTeX>
#endif

void
LALFindChirpBCVFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );

void
LALFindChirpBCVCFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );


#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPSPH_H */
