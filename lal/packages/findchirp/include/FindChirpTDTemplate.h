/*---------------------------------------------------------
 *      File Name: FindChirpTDTemplate.h
 *      author: S. Babak
 *      
 *----------------------------------------------------------
 */

#if 0

<lalVerbatim file="FindChirpTDTemplateHV">
Author: S. Babak
</lalVerbatim>

<lalLaTeX>
\section{Header \texttt{FindChirpTDTemplate.h}}
\label{s:FindChirpTDT.h}

This header file is a twin of FindChirpSP.h, FindChirpChisq.h and
. 
The main difference is that we filter data here using time domain waveforms.
User can choose any waveform from the inspiral package:
\\

TaylorT1, TaylorT2, TaylorT3, PadeT1, EOB

We try to use the same structures as in FindChirpSP.h

\begin{verbatim}
#include <lal/FindChirpSP.h>
\end{verbatim}

</lalLaTeX>
#endif



#ifndef _FINDCHIRPTDTEMPLATEH_H
#define _FINDCHIRPTDTEMPLATEH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpChisq.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPTDTH, "$Id$");

#define FINDCHIRPTDTEMPLATE_TLEN     1 
#define FINDCHIRPTDTEMPLATE_MSGETLEN "Template is too long for a given data\
					segment"

typedef struct
tagFindChirpChisqTDTParams{
   FindchirpChisqParams		*chisqParams,
   InspiralTemplate		*tmplt
}FindChirpChisqTDTParams;


void LALFindChirpTDTData(
   LALStatus 			*status,
   FindChirpSegmentVector	*fcSegVec,
   DataSegmentVector		*dataSegVec,
   FindChirpSPDataParams	*params
);

void LALFindChirpTDTVeto(
   LALStatus			*status,
   REAL4Vector			*chisqVec,
   FindChirpChisqInput		*input,
   FindChirpChisqTDTParams	*params
);

void LALFindChirpDstatVeto(
   LALStatus			*status,
   REAL4			*dist,
   FindChirpChisqInput		*input,
   FindChirpChisqParams		*params
);

void LALFindChirpTDTemplate(
   LALStatus			*status,
   FindChirpTemplate		*fcTmplt,
   InspiralTemplate		*tmplt,
   FindChirpSPDataParams	*params  /* note that params are different */
);

void LALFindChirpTDTFilterSegment(
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
);

#endif /* _FINDCHIRPTDTH_H */






