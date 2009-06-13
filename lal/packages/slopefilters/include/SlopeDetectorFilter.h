/*
*  Copyright (C) 2007 Edward Daw, Jolien Creighton
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

/*************************** <lalVerbatim file="SlopeDetectorFilterHV">
Author: Daw, E. J.
$Id$
**************************** </lalVerbatim> */

/*************************** <lalLaTeX>

\section{Header \texttt{SlopeDetectorFilter.h}}

Declares functions for slope detection. Defines the fields of
a structure \texttt{SLOPEFilterParams} whose fields are enumerated
and described below.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SlopeDetectorFilter.h>
\end{verbatim}

%[Generic documentation on the header; this is the main place to
%document any stuff not specific to the module]

\subsection*{Error conditions}
\input{SlopeDetectorFilterHE}

\subsection*{Structures}

\texttt{SLOPEFilterParams} is a structure defined in this header
which has the fields defined in table \ref{slopetable:paramstructure}.

\begin{table}
\begin{center}
\begin{tabular}{|l|c|r|} \hline
datatype & parameter & description \\ \hline\hline
\texttt{UINT4} & \texttt{forder} & number of bins $N$ used in filter \\ \hline
\texttt{REAL4*} & \texttt{tap} & filter taps \\ \hline
\texttt{UINT4*} & \texttt{history\_allocated} & 1 if set, 0 if not. \\ \hline
\texttt{UINT4*} & \texttt{taps\_set} & 1 if set, 0 if not. \\ \hline
\texttt{REAL4*} & \texttt{history} & history buffer \\ \hline
\texttt{UINT4} & \texttt{function\_select} & see function
descriptions \\ \hline
\texttt{REAL4} & \texttt{waveform\_offset} & offset of waveform
from bin, $0-1$ \\ \hline
\texttt{REAL4} & \texttt{sampling\_period\_s} & sampling period \\ \hline
\end{tabular}
\caption{Fields of the \texttt{SLOPEFilterParams} structure.}
\label{slopetable:paramstructure}
\end{center}
\end{table}

%[Document here any structures defined in the header.
%Also include any of them in the index; e.g.:]
% \index{\texttt{SlopeDetectorFilterOutput}}
% \index{\texttt{SlopeDetectorFilterInput}}
% \index{\texttt{SlopeDetectorFilterParams}}

\vfill{\footnotesize\input{SlopeDetectorFilterHV}}
\newpage\input{SlopeDetectorFilterC}
\newpage\input{SlopeDetectorFilterTestC}

********************************** </lalLaTeX> */

#ifndef _SLOPEDETECTORFILTER_H
#define _SLOPEDETECTORFILTER_H

#include <lal/LALStdlib.h>
/******* INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (SLOPEDETECTORFILTERH, "$Id$");

/******************************** <lalErrTable file="SlopeDetectorFilterHE"> */

#define SLOPEDETECTORFILTERH_EINPUTNULLP        1
#define SLOPEDETECTORFILTERH_EOUTPUTNULLP       2
#define SLOPEDETECTORFILTERH_ETAPSNULLP         3
#define SLOPEDETECTORFILTERH_EHISTNULLP         4
#define SLOPEDETECTORFILTERH_EINVFILTLEN        5
#define SLOPEDETECTORFILTERH_EDATATOOSHORT      6
#define SLOPEDETECTORFILTERH_EAMBHISTBIT        7
#define SLOPEDETECTORFILTERH_EINVALIDACTION     8
#define SLOPEDETECTORFILTERH_EBINOFFINVALID     9
#define SLOPEDETECTORFILTERH_EINVALIDTAPSBIT    10
#define SLOPEDETECTORFILTERH_EDIVBYZERO         11

#define SLOPEDETECTORFILTERH_MSGEINPUTNULLP   "Null input pointer"
#define SLOPEDETECTORFILTERH_MSGEOUTPUTNULLP  "Null output pointer"
#define SLOPEDETECTORFILTERH_MSGETAPSNULLP    "Inappropriate null taps pointer"
#define SLOPEDETECTORFILTERH_MSGEHISTNULLP    "Null history pointer"
#define SLOPEDETECTORFILTERH_MSGEINVFILTLEN   "Invalid filter length"
#define SLOPEDETECTORFILTERH_MSGEDATATOOSHORT "Data segment too short"
#define SLOPEDETECTORFILTERH_MSGEAMBHISTBIT "Ambiguous history buffer set bit"
#define SLOPEDETECTORFILTERH_MSGEINVALIDACTION "Filter action invalid"
#define SLOPEDETECTORFILTERH_MSGEBINOFFINVALID "Bin offset invalid"
#define SLOPEDETECTORFILTERH_MSGEINVALIDTAPSBIT "Invalid taps bit"
#define SLOPEDETECTORFILTERH_MSGEDIVBYZERO    "Division by zero"

/************************************ </lalErrTable> */

/****** DEFINE OTHER GLOBAL CONSTANTS OR MACROS ************/

#define MAX_SLOPE_DETECTOR_FILTER_ORDER    1000
#define GAUSSIAN_TAPS_NSIGMA_AT_EDGE       3.0
#define SLOPE_PI                           3.1415928
#define FILTER_OUTPUT_OFFSET               1
#define FILTER_OUTPUT_SLOPE                2
#define FILTER_OUTPUT_ALF                  3
#define FILTER_OUTPUT_TAPS_SET_BOXCAR      4
#define FILTER_OUTPUT_TAPS_SET_GAUSSIAN    5
#define FILTER_OUTPUT_TAPS_SET_SINE        6
#define FILTER_OUTPUT_TAPS_SET_USER        7
#define FILTER_OUTPUT_CONVOLVE             8

/****** DEFINE NEW STRUCTURES AND TYPES ************/

typedef struct tagSLOPEFilterParams{
  UINT4     forder;            /* number of taps in filter                */
  REAL4*    tap;               /* filter taps, may be null. Length ntaps  */
  UINT4*    history_allocated; /* 1 if yes, 0 if no                       */
  UINT4*    taps_set;          /* 1 if yes, 0 if no                       */
  REAL4*    history;           /* history buffer. Must be length ntaps-1  */
  UINT4     function_select;   /* 1=slope 2=offset 3=orsay magic combo    */
                               /* 4=setgaussian 5=setsinewave 6=setuser   */
                               /* 7=colvolve with existing taps           */
  REAL4     waveform_offset;   /* 0-1. Where does the center of the       */
                               /* waveform fall WRT the individual bins   */
  REAL4     sampling_period_s; /* sampling period in seconds              */
} SLOPEFilterParams;

/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

/****** FUNCTION DECLARATIONS ************************/

void
LALSlopeDetectorFilter( LALStatus            *status,
			REAL4Vector*         output_data,
			const REAL4Vector*   input_data,
			const UINT4          ntaps );

void
LALSlopeLineFitFilter( LALStatus                *status,
		       REAL4Vector*             output_data,
		       const REAL4Vector*       input_data,
		       const SLOPEFilterParams  fparams );

void
LALSlopeConvolutionFilter( LALStatus                *status,
			   REAL4Vector*             output_data,
			   const REAL4Vector*       input_data,
			   const SLOPEFilterParams  fparams );


#ifdef  __cplusplus
}
#endif

#endif /* _SLOPEDETECTORFILTER_H */






