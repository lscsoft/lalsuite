/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpTD.h
 *
 * Author: Brown, D. A., and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpTDHV">
Author: Brown, D. A., and Creighton, J. D. E.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{FindChirpTD.h}}
\label{s:FindChirpTD.h}

Provides structures and functions to condition interferometer data
and generate binary inspiral chirps using time domain waveforms.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpTD.h>
\end{verbatim}

</lalLaTeX>
#endif


#ifndef _FINDCHIRPTDH_H
#define _FINDCHIRPTDH_H

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


NRCSID (FINDCHIRPTDH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPTDH_ENULL 1
#define FINDCHIRPTDH_ENNUL 2
#define FINDCHIRPTDH_EALOC 3
#define FINDCHIRPTDH_ENUMZ 4
#define FINDCHIRPTDH_ESEGZ 5
#define FINDCHIRPTDH_EMISM 6
#define FINDCHIRPTDH_EDELT 7
#define FINDCHIRPTDH_EFLOW 8
#define FINDCHIRPTDH_EDYNR 9
#define FINDCHIRPTDH_EISTN 10
#define FINDCHIRPTDH_EDIVZ 11
#define FINDCHIRPTDH_EMAPX 12
#define FINDCHIRPTDH_ELONG 13
#define FINDCHIRPTDH_EEMTY 14
#define FINDCHIRPTDH_ESMPL 15
#define FINDCHIRPTDH_MSGENULL "Null pointer"
#define FINDCHIRPTDH_MSGENNUL "Non-null pointer"
#define FINDCHIRPTDH_MSGEALOC "Memory allocation error"
#define FINDCHIRPTDH_MSGENUMZ "Invalid number of segments"
#define FINDCHIRPTDH_MSGESEGZ "Invalid number of points in segments"
#define FINDCHIRPTDH_MSGEMISM "Mismatch between number of points in segments"
#define FINDCHIRPTDH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPTDH_MSGEFLOW "Low frequency cutoff is negative"
#define FINDCHIRPTDH_MSGEDYNR "Dynamic range scaling is zero or negative"
#define FINDCHIRPTDH_MSGEISTN "Truncation of inverse power spectrum is negative"
#define FINDCHIRPTDH_MSGEDIVZ "Attempting to divide by zero"
#define FINDCHIRPTDH_MSGEMAPX "Mismatch in waveform approximant"
#define FINDCHIRPTDH_MSGELONG "Time domain template too long"
#define FINDCHIRPTDH_MSGEEMTY "Could not find end of chirp in xfacVec"
#define FINDCHIRPTDH_MSGESMPL "Waveform sampling interval is too large"
/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Types}

None.
</lalLaTeX>
#endif


#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpTDHV}}
</lalLaTeX> 
#endif

#if 0
<lalLaTeX>
\newpage\input{FindChirpTDDataC}
</lalLaTeX>
#endif

void
LALFindChirpTDData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpTDTemplateC}
</lalLaTeX>
#endif

void
LALFindChirpTDTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *theTmplt,
    FindChirpTmpltParams       *params
    );

void
LALFindChirpTDNormalize(
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    FindChirpSegment           *fcSeg,
    FindChirpDataParams        *params
    );


#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPTDH_H */
