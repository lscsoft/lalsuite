/*---------------------------------------------------------
 *      File Name: FindChirpTDTemplate.h
 *      author: S. Babak
 *       
 *----------------------------------------------------------
 */

#if 0

<lalVerbatim file="FindChirpTDTemplateHV">
Author: S. Babak
$Id$ 
</lalVerbatim>

<lalLaTeX>
\section{Header \texttt{FindChirpTDTemplate.h}}
\label{s:FindChirpTDT.h}

This header file is a twin of FindChirpSP.h. 
It provides required structures and functions for filtering data
using time domain waveforms.
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

#define FINDCHIRPTDTEMPLATE_TLEN	1 
#define FINDCHIRPTDT_ETAPX 		2		
#define FINDCHIRPTDTEMPLATE_EFLOW	3
#define FINDCHIRPTDTEMPLATE_MSGETLEN "Template is too long for a given data segment"
#define FINDCHIRPTDT_MSGETAPX "Unknown approximant: must be TaylorT1, TaylorT3, PadeT1, EOB"
#define FINDCHIRPTDTEMPLATE_MSGEFLOW "Low frequency cutoff is different in fcTemplate and inspTemplate"

/*					
typedef
tagFindChirpChisqTDTParams{
   FindchirpChisqParams		*chisqParams,
   InspiralTemplate		*tmplt
}FindChirpChisqTDTParams;
*/

void LALFindChirpTDTData(
   LALStatus 			*status,
   FindChirpSegmentVector	*fcSegVec,
   DataSegmentVector		*dataSegVec,
   FindChirpDataParams		*params
);


void LALFindChirpTDTemplate(
   LALStatus			*status,
   FindChirpTemplate		*fcTmplt,
   InspiralTemplate		*tmplt,
   FindChirpDataParams		*params  /* note that params are different */
);

void LALFindChirpTDFilterSegment(
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
);
/*
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


*/
#endif /* _FINDCHIRPTDTH_H */






